-- Revert config_analysis_project_permission

BEGIN;

REVOKE ALL ON TABLE config.analysis_project FROM "gms-user";

COMMIT;
