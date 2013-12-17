-- Deploy disk_volume_permission
-- requires: disk_volume

BEGIN;

REVOKE ALL ON TABLE disk.volume FROM PUBLIC;
REVOKE ALL ON TABLE disk.volume FROM genome;
GRANT ALL ON TABLE disk.volume TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE disk.volume TO "gms-user";

COMMIT;
