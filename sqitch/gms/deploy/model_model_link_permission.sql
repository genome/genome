-- Deploy model_model_link_permission
-- requires: model_model_link

BEGIN;

REVOKE ALL ON TABLE model.model_link FROM PUBLIC;
REVOKE ALL ON TABLE model.model_link FROM genome;
GRANT ALL ON TABLE model.model_link TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.model_link TO "gms-user";

COMMIT;
