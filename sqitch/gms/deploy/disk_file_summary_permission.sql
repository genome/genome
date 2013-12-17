-- Deploy disk_file_summary_permission
-- requires: disk_file_summary

BEGIN;

REVOKE ALL ON TABLE disk.file_summary FROM PUBLIC;
REVOKE ALL ON TABLE disk.file_summary FROM genome;
GRANT ALL ON TABLE disk.file_summary TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE disk.file_summary TO "gms-user";

COMMIT;
