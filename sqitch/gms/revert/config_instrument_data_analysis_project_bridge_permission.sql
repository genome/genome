-- Revert config_instrument_data_analysis_project_bridge_permission

BEGIN;

REVOKE ALL ON TABLE config.instrument_data_analysis_project_bridge FROM "gms-user";

COMMIT;
