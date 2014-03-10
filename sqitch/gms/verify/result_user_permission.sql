-- Verify result_user_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'result."user"', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'result."user"', 'SELECT')::int;

ROLLBACK;
