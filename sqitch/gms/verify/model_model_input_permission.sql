-- Verify model_model_input_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.model_input', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.model_input', 'SELECT')::int;

ROLLBACK;
