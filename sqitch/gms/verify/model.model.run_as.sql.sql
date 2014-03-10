-- Verify model.model.run_as.sql

BEGIN;

    SELECT run_as FROM model.model WHERE FALSE;

ROLLBACK;
