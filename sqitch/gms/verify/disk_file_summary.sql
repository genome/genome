-- Verify disk_file_summary

BEGIN;

SELECT id, updated_at, created_at, created_by, allocation_id, file,
    digest, size_in_bytes, is_symlink, destination
FROM disk.file_summary
WHERE FALSE;

ROLLBACK;
