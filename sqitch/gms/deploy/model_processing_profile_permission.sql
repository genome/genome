-- Deploy model_processing_profile_permission
-- requires: model_processing_profile

BEGIN;

REVOKE ALL ON TABLE model.processing_profile FROM PUBLIC;
REVOKE ALL ON TABLE model.processing_profile FROM genome;
GRANT ALL ON TABLE model.processing_profile TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.processing_profile TO "gms-user";

COMMIT;
