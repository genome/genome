-- Verify workflow_schema

BEGIN;

SELECT pg_catalog.has_schema_privilege('genome', 'workflow', 'usage');

ROLLBACK;
