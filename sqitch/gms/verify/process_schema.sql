-- Verify process_schema

BEGIN;

SELECT pg_catalog.has_schema_privilege('process', 'usage');

ROLLBACK;
