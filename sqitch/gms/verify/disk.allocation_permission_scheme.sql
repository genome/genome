-- Verify disk.allocation_permission_scheme

BEGIN;

SELECT id FROM disk.allocation_permission_scheme WHERE FALSE;
SELECT permission_scheme_id FROM disk.allocation WHERE FALSE;

ROLLBACK;
