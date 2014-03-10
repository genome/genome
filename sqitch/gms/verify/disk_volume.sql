-- Verify disk_volume

BEGIN;

SELECT id, hostname,physical_path,mount_path, total_kb,
    unallocated_kb, disk_status, can_allocate, doubles_space
FROM disk.volume
WHERE FALSE;

ROLLBACK;
