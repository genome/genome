-- Verify instrument_data_attribute_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'instrument.data_attribute', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'instrument.data_attribute', 'SELECT')::int;

ROLLBACK;
