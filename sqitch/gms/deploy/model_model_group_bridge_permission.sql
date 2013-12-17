-- Deploy model_model_group_bridge_permission
-- requires: model_model_group_bridge

BEGIN;

REVOKE ALL ON TABLE model.model_group_bridge FROM PUBLIC;
REVOKE ALL ON TABLE model.model_group_bridge FROM genome;
GRANT ALL ON TABLE model.model_group_bridge TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.model_group_bridge TO "gms-user";

COMMIT;
