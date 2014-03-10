-- Deploy disk_volume_group_bridge_permission
-- requires: disk_volume_group_bridge

BEGIN;

REVOKE ALL ON TABLE disk.volume_group_bridge FROM PUBLIC;
REVOKE ALL ON TABLE disk.volume_group_bridge FROM genome;
GRANT ALL ON TABLE disk.volume_group_bridge TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE disk.volume_group_bridge TO "gms-user";

COMMIT;
