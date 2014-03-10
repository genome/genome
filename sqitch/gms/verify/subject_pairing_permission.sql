-- Verify subject_pairing_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'subject.pairing', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'subject.pairing', 'SELECT')::int;

ROLLBACK;
