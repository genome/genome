-- Verify model.model.created_by

BEGIN;

SELECT created_by FROM model.model WHERE FALSE;

ROLLBACK;
