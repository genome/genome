-- Verify add_config_id_to_anp_model_bridges

BEGIN;
SELECT profile_item_id FROM config.analysis_project_model_bridge WHERE false;
ROLLBACK;
