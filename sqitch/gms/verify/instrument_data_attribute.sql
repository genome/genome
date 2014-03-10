-- Verify instrument_data_attribute

BEGIN;

SELECT instrument_data_id, attribute_label, attribute_value, nomenclature
FROM instrument.data_attribute
WHERE FALSE;

ROLLBACK;
