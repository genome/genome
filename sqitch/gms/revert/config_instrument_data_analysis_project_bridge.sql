-- Revert config_instrument_data_analysis_project_bridge

BEGIN;

DROP TABLE config.instrument_data_analysis_project_bridge;

COMMIT;
