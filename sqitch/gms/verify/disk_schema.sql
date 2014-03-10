-- Verify config_schema

BEGIN;

SELECT pg_catalog.has_schema_privilege('disk', 'usage');

ROLLBACK;
