-- Deploy model_processing_profile_param_permission
-- requires: model_processing_profile_param

BEGIN;

REVOKE ALL ON TABLE model.processing_profile_param FROM PUBLIC;
REVOKE ALL ON TABLE model.processing_profile_param FROM genome;
GRANT ALL ON TABLE model.processing_profile_param TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.processing_profile_param TO "gms-user";

COMMIT;
