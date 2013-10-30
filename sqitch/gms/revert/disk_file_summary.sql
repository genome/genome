-- Revert disk_file_summary

BEGIN;

DROP TABLE IF EXISTS disk.file_summary;

COMMIT;
