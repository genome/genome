-- Revert disk_file_summary_permission

BEGIN;

REVOKE ALL ON TABLE disk.file_summary FROM "gms-user";

COMMIT;
