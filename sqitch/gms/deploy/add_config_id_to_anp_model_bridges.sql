-- Deploy add_config_id_to_anp_model_bridges
-- requires: config_analysis_project_model_bridge

BEGIN;
ALTER TABLE config.analysis_project_model_bridge ADD COLUMN profile_item_id character varying(64) REFERENCES config.profile_item(id);
COMMIT;
