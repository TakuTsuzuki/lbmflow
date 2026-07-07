use serde_json::{json, Value};

pub fn bioprocess_schema_json() -> String {
    serde_json::to_string_pretty(&bioprocess_schema()).expect("schema JSON must serialize")
}

fn bioprocess_schema() -> Value {
    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "BioprocessScenario",
        "type": "object",
        "additionalProperties": false,
        "required": [
            "version",
            "name",
            "credibility_tier",
            "reactor",
            "fluids",
            "operation",
            "physics",
            "qoi",
            "run",
            "outputs"
        ],
        "properties": {
            "version": { "type": "string", "const": "bioprocess-1.0" },
            "name": { "type": "string" },
            "credibility_tier": {
                "type": "string",
                "enum": ["screening", "engineering", "evidence"]
            },
            "reactor": { "$ref": "#/$defs/ReactorSpec" },
            "fluids": { "$ref": "#/$defs/FluidsSpec" },
            "operation": { "$ref": "#/$defs/OperationSpec" },
            "physics": {
                "oneOf": [
                    { "$ref": "#/$defs/PhysicsModel" },
                    {
                        "type": "array",
                        "items": { "$ref": "#/$defs/PhysicsModel" },
                        "minItems": 1
                    }
                ]
            },
            "cells": { "$ref": "#/$defs/CellsSpec" },
            "qoi": { "$ref": "#/$defs/QoiSpec" },
            "run": { "$ref": "#/$defs/RunSpec" },
            "outputs": { "$ref": "#/$defs/OutputsSpec" }
        },
        "$defs": {
            "ReactorSpec": {
                "type": "object",
                "additionalProperties": false,
                "required": [
                    "kind",
                    "vessel_diameter_m",
                    "liquid_height_m",
                    "working_volume_m3",
                    "impellers",
                    "baffles",
                    "spargers"
                ],
                "properties": {
                    "kind": { "const": "stirred_tank" },
                    "vessel_diameter_m": { "type": "number" },
                    "liquid_height_m": { "type": "number" },
                    "working_volume_m3": { "type": "number" },
                    "impellers": {
                        "type": "array",
                        "items": { "$ref": "#/$defs/ImpellerSpec" }
                    },
                    "baffles": {
                        "type": "array",
                        "items": { "$ref": "#/$defs/BaffleSpec" }
                    },
                    "spargers": {
                        "type": "array",
                        "items": { "$ref": "#/$defs/SpargerSpec" }
                    }
                }
            },
            "ImpellerSpec": {
                "type": "object",
                "additionalProperties": false,
                "required": [
                    "kind",
                    "diameter_m",
                    "clearance_from_bottom_m",
                    "rotational_speed_rpm",
                    "blade_count"
                ],
                "properties": {
                    "kind": {
                        "type": "string",
                        "enum": ["rushton", "pitched_blade", "marine", "custom_marker_set"]
                    },
                    "diameter_m": { "type": "number" },
                    "clearance_from_bottom_m": { "type": "number" },
                    "rotational_speed_rpm": { "type": "number" },
                    "blade_count": { "type": "integer", "minimum": 0 }
                }
            },
            "BaffleSpec": {
                "type": "object",
                "additionalProperties": false,
                "required": ["count", "width_m", "thickness_m", "wall_attached"],
                "properties": {
                    "count": { "type": "integer", "minimum": 0 },
                    "width_m": { "type": "number" },
                    "thickness_m": { "type": "number" },
                    "wall_attached": { "type": "boolean" }
                }
            },
            "SpargerSpec": {
                "oneOf": [
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": [
                            "kind",
                            "center_z_m",
                            "outer_radius_m",
                            "orifice_count",
                            "orifice_diameter_m",
                            "gas_volumetric_flow_m3_per_s",
                            "vvm",
                            "inlet_phase"
                        ],
                        "properties": {
                            "kind": { "const": "ring" },
                            "center_z_m": { "type": "number" },
                            "outer_radius_m": { "type": "number" },
                            "orifice_count": { "type": "integer", "minimum": 0 },
                            "orifice_diameter_m": { "type": "number" },
                            "gas_volumetric_flow_m3_per_s": { "type": ["number", "null"] },
                            "vvm": { "type": ["number", "null"] },
                            "inlet_phase": { "type": "string" }
                        }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": [
                            "kind",
                            "center_z_m",
                            "length_m",
                            "diameter_m",
                            "gas_volumetric_flow_m3_per_s",
                            "vvm",
                            "inlet_phase"
                        ],
                        "properties": {
                            "kind": { "const": "pipe" },
                            "center_z_m": { "type": "number" },
                            "length_m": { "type": "number" },
                            "diameter_m": { "type": "number" },
                            "gas_volumetric_flow_m3_per_s": { "type": ["number", "null"] },
                            "vvm": { "type": ["number", "null"] },
                            "inlet_phase": { "type": "string" }
                        }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": [
                            "kind",
                            "center_z_m",
                            "positions",
                            "gas_volumetric_flow_m3_per_s",
                            "vvm",
                            "inlet_phase"
                        ],
                        "properties": {
                            "kind": { "const": "point_orifices" },
                            "center_z_m": { "type": "number" },
                            "positions": {
                                "type": "array",
                                "items": {
                                    "type": "array",
                                    "prefixItems": [
                                        { "type": "number" },
                                        { "type": "number" },
                                        { "type": "number" }
                                    ],
                                    "items": false
                                }
                            },
                            "gas_volumetric_flow_m3_per_s": { "type": ["number", "null"] },
                            "vvm": { "type": ["number", "null"] },
                            "inlet_phase": { "type": "string" }
                        }
                    }
                ]
            },
            "FluidsSpec": {
                "type": "object",
                "additionalProperties": false,
                "required": [
                    "liquid_density_kg_m3",
                    "liquid_viscosity_pa_s",
                    "gas_density_kg_m3",
                    "gas_viscosity_pa_s",
                    "surface_tension_n_m",
                    "oxygen_diffusivity_m2_per_s",
                    "henry_constant"
                ],
                "properties": {
                    "liquid_density_kg_m3": { "type": "number" },
                    "liquid_viscosity_pa_s": { "type": "number" },
                    "gas_density_kg_m3": { "type": ["number", "null"] },
                    "gas_viscosity_pa_s": { "type": ["number", "null"] },
                    "surface_tension_n_m": { "type": ["number", "null"] },
                    "oxygen_diffusivity_m2_per_s": { "type": ["number", "null"] },
                    "henry_constant": { "type": ["number", "null"] }
                }
            },
            "OperationSpec": {
                "type": "object",
                "additionalProperties": false,
                "required": ["duration_s", "gas_inlet_temp_c", "initial_condition"],
                "properties": {
                    "duration_s": { "type": "number" },
                    "gas_inlet_temp_c": { "type": ["number", "null"] },
                    "initial_condition": {
                        "oneOf": [
                            {
                                "type": "object",
                                "additionalProperties": false,
                                "required": ["kind"],
                                "properties": { "kind": { "const": "quiescent" } }
                            },
                            {
                                "type": "object",
                                "additionalProperties": false,
                                "required": ["kind", "path"],
                                "properties": {
                                    "kind": { "const": "existing_checkpoint" },
                                    "path": { "type": "string" }
                                }
                            }
                        ]
                    }
                }
            },
            "PhysicsModel": {
                "oneOf": [
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind"],
                        "properties": { "kind": { "const": "single_phase" } }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "interface_width_m", "mobility_m2_per_s", "contact_angle_deg"],
                        "properties": {
                            "kind": { "const": "resolved_phase_field" },
                            "interface_width_m": { "type": "number" },
                            "mobility_m2_per_s": { "type": "number" },
                            "contact_angle_deg": { "type": ["number", "null"] }
                        }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "max_bubble_count"],
                        "properties": {
                            "kind": { "const": "point_bubble" },
                            "max_bubble_count": { "type": "integer", "minimum": 0 }
                        }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "diffusivity_m2_per_s", "initial_pulse"],
                        "properties": {
                            "kind": { "const": "passive_scalar" },
                            "diffusivity_m2_per_s": { "type": "number" },
                            "initial_pulse": { "type": ["object", "null"] }
                        }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "henry_constant", "interfacial_flux_model", "our_model"],
                        "properties": {
                            "kind": { "const": "oxygen" },
                            "henry_constant": { "type": "number" },
                            "interfacial_flux_model": {
                                "type": "string",
                                "enum": ["henry_equilibrium", "constant_kl", "interfacial_area"]
                            },
                            "our_model": {
                                "type": "string",
                                "enum": ["none", "constant", "monod"]
                            }
                        }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": [
                            "kind",
                            "count",
                            "seed",
                            "record_shear",
                            "record_oxygen",
                            "damage_model"
                        ],
                        "properties": {
                            "kind": { "const": "cell_tracer" },
                            "count": { "type": "integer", "minimum": 0 },
                            "seed": { "type": "integer", "minimum": 0 },
                            "record_shear": { "type": "boolean" },
                            "record_oxygen": { "type": "boolean" },
                            "damage_model": { "type": ["object", "null"] }
                        }
                    }
                ]
            },
            "CellsSpec": {
                "type": "object",
                "additionalProperties": false,
                "required": ["count", "seed", "record_shear", "record_oxygen", "damage_model"],
                "properties": {
                    "count": { "type": "integer", "minimum": 0 },
                    "seed": { "type": "integer", "minimum": 0 },
                    "record_shear": { "type": "boolean" },
                    "record_oxygen": { "type": "boolean" },
                    "damage_model": { "type": ["object", "null"] }
                }
            },
            "QoiSpec": {
                "type": "object",
                "additionalProperties": false,
                "required": [
                    "power",
                    "mixing",
                    "gas_holdup",
                    "bubble_size",
                    "kla",
                    "shear_exposure",
                    "oxygen_exposure",
                    "calibration_dataset_id",
                    "holdout_dataset_id"
                ],
                "properties": {
                    "power": { "type": ["object", "null"], "additionalProperties": false },
                    "mixing": { "type": ["object", "null"], "additionalProperties": false },
                    "gas_holdup": { "type": ["object", "null"], "additionalProperties": false },
                    "bubble_size": { "type": ["object", "null"], "additionalProperties": false },
                    "kla": { "type": ["object", "null"], "additionalProperties": false },
                    "shear_exposure": { "type": ["object", "null"], "additionalProperties": false },
                    "oxygen_exposure": { "type": ["object", "null"], "additionalProperties": false },
                    "calibration_dataset_id": { "type": ["string", "null"] },
                    "holdout_dataset_id": { "type": ["string", "null"] }
                }
            },
            "RunSpec": {
                "type": "object",
                "additionalProperties": false,
                "required": ["steps", "dt_s", "grid_nx", "grid_ny", "grid_nz", "backend", "precision", "lattice"],
                "properties": {
                    "steps": { "type": "integer", "minimum": 0 },
                    "dt_s": { "type": "number" },
                    "grid_nx": { "type": "integer", "minimum": 1 },
                    "grid_ny": { "type": "integer", "minimum": 1 },
                    "grid_nz": { "type": "integer", "minimum": 1 },
                    "backend": { "type": ["string", "null"], "enum": ["auto", "cpu", "gpu", null] },
                    "precision": { "type": ["string", "null"], "enum": ["f32", "f64", null] },
                    "lattice": { "type": ["string", "null"], "enum": ["d2q9", "d3q19", "d3q27", null] }
                }
            },
            "UnitReport": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "lattice": { "type": "object" },
                    "groups": { "type": "object" },
                    "feasibility": { "type": "object" },
                    "matching_priority": { "type": "object" }
                }
            },
            "OutputsSpec": {
                "type": "object",
                "additionalProperties": false,
                "required": [
                    "manifest_path",
                    "fields_every_n_steps",
                    "probes_every_n_steps",
                    "emit_qoi_json",
                    "emit_qoi_csv"
                ],
                "properties": {
                    "manifest_path": { "type": "string" },
                    "fields_every_n_steps": { "type": ["integer", "null"], "minimum": 0 },
                    "probes_every_n_steps": { "type": ["integer", "null"], "minimum": 0 },
                    "emit_qoi_json": { "type": "boolean" },
                    "emit_qoi_csv": { "type": "boolean" }
                }
            }
        }
    })
}
