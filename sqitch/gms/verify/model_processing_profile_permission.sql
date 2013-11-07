-- Verify model_processing_profile_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.processing_profile', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.processing_profile', 'SELECT')::int;

ROLLBACK;
