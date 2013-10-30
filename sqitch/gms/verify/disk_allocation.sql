-- Verify disk_allocation

BEGIN;

SELECT id, allocation_path, disk_group_name, group_subdirectory,
    mount_path, kilobytes_requested, kilobytes_used, owner_class_name,
    owner_id, creation_time, reallocation_time, original_kilobytes_requested,
    archivable, preserved, kilobytes_used_time, archive_after_time, status
FROM disk.allocation
WHERE FALSE;

ROLLBACK;
