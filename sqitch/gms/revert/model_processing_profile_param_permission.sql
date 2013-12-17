-- Revert model_processing_profile_param_permission

BEGIN;

REVOKE ALL ON TABLE model.processing_profile_param FROM "gms-user";

COMMIT;
