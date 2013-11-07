-- Verify instrument_fragment_library_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'instrument.fragment_library', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'instrument.fragment_library', 'SELECT')::int;

ROLLBACK;
