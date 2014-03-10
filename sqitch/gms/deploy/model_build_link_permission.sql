-- Deploy model_build_link_permission
-- requires: model_build_link

BEGIN;

REVOKE ALL ON TABLE model.build_link FROM PUBLIC;
REVOKE ALL ON TABLE model.build_link FROM genome;
GRANT ALL ON TABLE model.build_link TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.build_link TO "gms-user";

COMMIT;
