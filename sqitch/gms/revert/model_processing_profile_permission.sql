-- Revert model_processing_profile_permission

BEGIN;

REVOKE ALL ON TABLE model.processing_profile FROM "gms-user";

COMMIT;
