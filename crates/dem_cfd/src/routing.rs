//! Spatial route construction owned by the DEM-CFD coupling package.
//!
//! FIELD supplies partition ownership, SOIL supplies particle state and DEM
//! ownership, and this module turns those facts into opaque GRASS routes.

use field_core::PartitionDirectory;
use grass_multi::{EntityId, ReceivedPayload, RoutedPayload};
use std::collections::BTreeMap;
use std::fmt;

const PARTICLE_BYTES: usize = 8 + 4 + 3 * 8 + 3 * 8 + 8 + 3 * 8;
const FORCE_BYTES: usize = 8 + 4 + 3 * 8;

/// Particle state required by the unresolved DEM-CFD spatial seam.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RoutedParticle {
    /// Stable particle ID across migration and repartitioning.
    pub id: u64,
    /// Current owner within the DEM role communicator.
    pub dem_owner: i32,
    /// Particle center in world coordinates.
    pub center: [f64; 3],
    /// Particle velocity in world coordinates.
    pub velocity: [f64; 3],
    /// Physical particle radius.
    pub radius: f64,
    /// Fluid force the particle carried out of the *previous* coupling step,
    /// travelling with the particle across CFD-partition and DEM-ownership
    /// boundaries. The receiving CFD owner needs it to form the trapezoidal
    /// (velocity-Verlet) momentum-sink reaction `0.5·(F_prev + F_now)` even on
    /// the first step it sees a freshly migrated particle; otherwise a boundary
    /// crossing would silently degrade the sink to a rectangle rule and break
    /// exact cross-role impulse conservation. Zero on the primed first step.
    pub prev_force: [f64; 3],
}

/// Force result returned by a CFD owner to the particle's DEM owner.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RoutedForce {
    /// Stable particle ID.
    pub id: u64,
    /// DEM role-rank that owns the particle.
    pub dem_owner: i32,
    /// Total fluid force on the particle.
    pub force: [f64; 3],
}

/// Build containing-partition routes for particle kinematics.
pub fn route_particles(
    directory: &PartitionDirectory,
    particles: &[RoutedParticle],
) -> Result<Vec<RoutedPayload>, ParticleRoutingError> {
    particles
        .iter()
        .map(|particle| {
            let destination = directory.owner_rank(particle.center).ok_or(
                ParticleRoutingError::OutsideFieldDomain {
                    id: particle.id,
                    center: particle.center,
                },
            )?;
            Ok(RoutedPayload::new(
                destination,
                EntityId(particle.id),
                encode_particle(*particle),
            ))
        })
        .collect()
}

/// Decode particle records delivered by GRASS, verifying frame/entity IDs.
pub fn decode_particles(
    incoming: &[ReceivedPayload],
) -> Result<Vec<RoutedParticle>, ParticleRoutingError> {
    incoming
        .iter()
        .map(|record| {
            let particle = decode_particle(&record.payload)?;
            if particle.id != record.entity_id.0 {
                return Err(ParticleRoutingError::EntityIdMismatch {
                    envelope: record.entity_id.0,
                    payload: particle.id,
                });
            }
            Ok(particle)
        })
        .collect()
}

/// Route CFD force results back to current DEM owners.
pub fn route_forces(forces: &[RoutedForce]) -> Vec<RoutedPayload> {
    forces
        .iter()
        .map(|force| RoutedPayload::new(force.dem_owner, EntityId(force.id), encode_force(*force)))
        .collect()
}

/// Decode and deterministically reduce force contributions by particle ID.
/// Multiple CFD owners may contribute when a future finite-support mapping
/// routes one particle to several partitions.
pub fn reduce_forces(
    incoming: &[ReceivedPayload],
) -> Result<Vec<RoutedForce>, ParticleRoutingError> {
    let mut reduced: BTreeMap<(u64, i32), [f64; 3]> = BTreeMap::new();
    for record in incoming {
        let force = decode_force(&record.payload)?;
        if force.id != record.entity_id.0 {
            return Err(ParticleRoutingError::EntityIdMismatch {
                envelope: record.entity_id.0,
                payload: force.id,
            });
        }
        let total = reduced.entry((force.id, force.dem_owner)).or_default();
        for (total_axis, contribution) in total.iter_mut().zip(force.force) {
            *total_axis += contribution;
        }
    }
    Ok(reduced
        .into_iter()
        .map(|((id, dem_owner), force)| RoutedForce {
            id,
            dem_owner,
            force,
        })
        .collect())
}

fn encode_particle(particle: RoutedParticle) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(PARTICLE_BYTES);
    bytes.extend_from_slice(&particle.id.to_le_bytes());
    bytes.extend_from_slice(&particle.dem_owner.to_le_bytes());
    for value in particle.center.into_iter().chain(particle.velocity) {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    bytes.extend_from_slice(&particle.radius.to_le_bytes());
    for value in particle.prev_force {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    bytes
}

fn decode_particle(bytes: &[u8]) -> Result<RoutedParticle, ParticleRoutingError> {
    if bytes.len() != PARTICLE_BYTES {
        return Err(ParticleRoutingError::MalformedPayload {
            expected: PARTICLE_BYTES,
            actual: bytes.len(),
        });
    }
    let mut cursor = 0;
    let id = read_u64(bytes, &mut cursor);
    let dem_owner = read_i32(bytes, &mut cursor);
    let center = std::array::from_fn(|_| read_f64(bytes, &mut cursor));
    let velocity = std::array::from_fn(|_| read_f64(bytes, &mut cursor));
    let radius = read_f64(bytes, &mut cursor);
    let prev_force = std::array::from_fn(|_| read_f64(bytes, &mut cursor));
    Ok(RoutedParticle {
        id,
        dem_owner,
        center,
        velocity,
        radius,
        prev_force,
    })
}

fn encode_force(force: RoutedForce) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(FORCE_BYTES);
    bytes.extend_from_slice(&force.id.to_le_bytes());
    bytes.extend_from_slice(&force.dem_owner.to_le_bytes());
    for value in force.force {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    bytes
}

fn decode_force(bytes: &[u8]) -> Result<RoutedForce, ParticleRoutingError> {
    if bytes.len() != FORCE_BYTES {
        return Err(ParticleRoutingError::MalformedPayload {
            expected: FORCE_BYTES,
            actual: bytes.len(),
        });
    }
    let mut cursor = 0;
    let id = read_u64(bytes, &mut cursor);
    let dem_owner = read_i32(bytes, &mut cursor);
    let force = std::array::from_fn(|_| read_f64(bytes, &mut cursor));
    Ok(RoutedForce {
        id,
        dem_owner,
        force,
    })
}

fn read_u64(bytes: &[u8], cursor: &mut usize) -> u64 {
    let value = u64::from_le_bytes(bytes[*cursor..*cursor + 8].try_into().unwrap());
    *cursor += 8;
    value
}

fn read_i32(bytes: &[u8], cursor: &mut usize) -> i32 {
    let value = i32::from_le_bytes(bytes[*cursor..*cursor + 4].try_into().unwrap());
    *cursor += 4;
    value
}

fn read_f64(bytes: &[u8], cursor: &mut usize) -> f64 {
    let value = f64::from_le_bytes(bytes[*cursor..*cursor + 8].try_into().unwrap());
    *cursor += 8;
    value
}

/// Invalid coupling-owned particle route or payload.
#[derive(Debug, Clone, PartialEq)]
pub enum ParticleRoutingError {
    /// Particle center has no FIELD owner.
    OutsideFieldDomain {
        /// Particle ID.
        id: u64,
        /// Rejected center.
        center: [f64; 3],
    },
    /// Envelope and payload did not describe the same entity.
    EntityIdMismatch {
        /// GRASS envelope ID.
        envelope: u64,
        /// Coupling payload ID.
        payload: u64,
    },
    /// Fixed-size coupling payload was truncated or extended.
    MalformedPayload {
        /// Required byte length.
        expected: usize,
        /// Received byte length.
        actual: usize,
    },
}

impl fmt::Display for ParticleRoutingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::OutsideFieldDomain { id, center } => {
                write!(f, "particle {id} at {center:?} has no FIELD owner")
            }
            Self::EntityIdMismatch { envelope, payload } => write!(
                f,
                "routed entity ID mismatch: envelope {envelope}, payload {payload}"
            ),
            Self::MalformedPayload { expected, actual } => write!(
                f,
                "malformed routed DEM-CFD payload: expected {expected} bytes, received {actual}"
            ),
        }
    }
}

impl std::error::Error for ParticleRoutingError {}

#[cfg(test)]
mod tests {
    use super::*;
    use field_core::{UniformMeshConfig, Vec3};

    fn directory() -> PartitionDirectory {
        PartitionDirectory::from_uniform_config(
            &UniformMeshConfig {
                nx: 10,
                ny: 1,
                nz: 1,
                ng: 1,
                bounds_lo: [0.0; 3],
                bounds_hi: [1.0; 3],
                y_edges: None,
                z_edges: None,
            },
            [2, 1, 1],
        )
    }

    fn particle(id: u64, owner: i32, x: f64) -> RoutedParticle {
        RoutedParticle {
            id,
            dem_owner: owner,
            center: [x, 0.5, 0.5],
            velocity: [x, -x, 0.0],
            radius: 0.01,
            prev_force: [0.5 * x, -0.25 * x, 0.125],
        }
    }

    #[test]
    fn routes_particle_centers_to_field_owners() {
        let routes =
            route_particles(&directory(), &[particle(1, 0, 0.1), particle(2, 1, 0.5)]).unwrap();
        assert_eq!(routes[0].destination, 0);
        assert_eq!(routes[1].destination, 1);
    }

    #[test]
    fn particle_and_force_codecs_round_trip_and_reduce() {
        let p = particle(8, 2, 0.75);
        assert_eq!(decode_particle(&encode_particle(p)).unwrap(), p);
        let incoming = vec![
            ReceivedPayload {
                source: 0,
                entity_id: EntityId(8),
                payload: encode_force(RoutedForce {
                    id: 8,
                    dem_owner: 2,
                    force: [1.0, 2.0, 3.0],
                }),
            },
            ReceivedPayload {
                source: 1,
                entity_id: EntityId(8),
                payload: encode_force(RoutedForce {
                    id: 8,
                    dem_owner: 2,
                    force: [4.0, 5.0, 6.0],
                }),
            },
        ];
        assert_eq!(
            reduce_forces(&incoming).unwrap(),
            vec![RoutedForce {
                id: 8,
                dem_owner: 2,
                force: [5.0, 7.0, 9.0],
            }]
        );
    }

    #[test]
    fn outside_particle_fails_closed() {
        let error = route_particles(&directory(), &[particle(4, 0, 1.1)]).unwrap_err();
        assert!(matches!(
            error,
            ParticleRoutingError::OutsideFieldDomain { id: 4, .. }
        ));
    }

    #[allow(dead_code)]
    fn _vec3_type_is_shared(_: Vec3) {}
}
