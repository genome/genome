-- Verify disk.allocation_permission_scheme_name

BEGIN;

SELECT name FROM disk.allocation_permission_scheme WHERE FALSE;

ROLLBACK;
