-- Verify instrument.fragment_library_unique_name

BEGIN;

SELECT 1/COUNT(*) FROM pg_constraint AS c
INNER JOIN pg_namespace AS n ON c.connamespace = n.oid
WHERE n.nspname = 'instrument' AND c.conname = 'fragment_library_unique_full_name';

ROLLBACK;
