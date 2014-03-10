-- Deploy instrument_data_attribute
-- requires: instrument_data

BEGIN;

CREATE TABLE IF NOT EXISTS instrument.data_attribute (
    instrument_data_id character varying(64) NOT NULL,
    attribute_label character varying(64) NOT NULL,
    attribute_value character varying(512) NOT NULL,
    nomenclature character varying(64) NOT NULL,
    CONSTRAINT data_attribute_pkey PRIMARY KEY (instrument_data_id, attribute_label, attribute_value, nomenclature),
    CONSTRAINT data_attribute_instrument_data_id_fkey FOREIGN KEY (instrument_data_id) REFERENCES instrument.data(id)
);

COMMIT;
