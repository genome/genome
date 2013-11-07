-- Deploy model_model_group_permission
-- requires: model_model_group

BEGIN;

REVOKE ALL ON TABLE model.model_group FROM PUBLIC;
REVOKE ALL ON TABLE model.model_group FROM genome;
GRANT ALL ON TABLE model.model_group TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.model_group TO "gms-user";

COMMIT;
