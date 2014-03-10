-- Deploy config_instrument_data_analysis_project_bridge_permission
-- requires: config_instrument_data_analysis_project_bridge

BEGIN;

REVOKE ALL ON TABLE config.instrument_data_analysis_project_bridge FROM PUBLIC;
REVOKE ALL ON TABLE config.instrument_data_analysis_project_bridge FROM genome;
GRANT ALL ON TABLE config.instrument_data_analysis_project_bridge TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.instrument_data_analysis_project_bridge TO "gms-user";

COMMIT;
