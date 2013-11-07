-- Verify instrument_data_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'instrument.data', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'instrument.data', 'SELECT')::int;

ROLLBACK;
