-- Verify timeline_allocation

BEGIN;

SELECT kilobytes_requested, absolute_path, status
FROM timeline.allocation
WHERE FALSE;

ROLLBACK;
