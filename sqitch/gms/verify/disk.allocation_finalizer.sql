-- Verify disk.allocation_finalizer

BEGIN;

SELECT id FROM disk.allocation_finalizer WHERE FALSE;
SELECT finalizer_id FROM disk.allocation WHERE FALSE;

ROLLBACK;
