-- Deploy model_model_permission
-- requires: model_model

BEGIN;

REVOKE ALL ON TABLE model.model FROM PUBLIC;
REVOKE ALL ON TABLE model.model FROM genome;
GRANT ALL ON TABLE model.model TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.model TO "gms-user";

COMMIT;
