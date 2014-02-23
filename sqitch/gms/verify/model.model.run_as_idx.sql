-- Verify model.model.run_as_idx

BEGIN;

SELECT 1/count(*) FROM pg_class AS c
    INNER JOIN pg_namespace AS n ON c.relnamespace = n.oid
    WHERE n.nspname = 'model' AND c.relkind = 'i' and c.relname = 'model__run_as_idx';

ROLLBACK;
