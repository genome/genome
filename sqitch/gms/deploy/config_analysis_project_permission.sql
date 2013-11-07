-- Deploy config_analysis_project_permission
-- requires: config_analysis_project

BEGIN;

REVOKE ALL ON TABLE config.analysis_project FROM PUBLIC;
REVOKE ALL ON TABLE config.analysis_project FROM genome;
GRANT ALL ON TABLE config.analysis_project TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.analysis_project TO "gms-user";
COMMIT;
