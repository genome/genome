-- Revert result.input.index_input_value_input_name

BEGIN;

DROP INDEX result.sri_iniv;

COMMIT;
