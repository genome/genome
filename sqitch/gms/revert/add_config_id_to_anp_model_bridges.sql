-- Revert add_config_id_to_anp_model_bridges

BEGIN;
ALTER TABLE config.analysis_project_model_bridge DROP COLUMN profile_item_id;
COMMIT;
