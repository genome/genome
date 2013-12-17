-- Deploy disk_allocation_permission
-- requires: disk_allocation

BEGIN;

REVOKE ALL ON TABLE disk.allocation FROM PUBLIC;
REVOKE ALL ON TABLE disk.allocation FROM genome;
GRANT ALL ON TABLE disk.allocation TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE disk.allocation TO "gms-user";

COMMIT;
