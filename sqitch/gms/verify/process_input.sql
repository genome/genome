-- Verify process_input

BEGIN;

SELECT id, process_id, value_class_name, value_id, name, array_index
FROM process.input WHERE FALSE;

ROLLBACK;
