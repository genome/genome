-- Deploy timeline_base_permission
-- requires: timeline_base

BEGIN;

REVOKE ALL ON TABLE timeline.base FROM PUBLIC;
REVOKE ALL ON TABLE timeline.base FROM genome;
GRANT ALL ON TABLE timeline.base TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE timeline.base TO "gms-user";

COMMIT;
