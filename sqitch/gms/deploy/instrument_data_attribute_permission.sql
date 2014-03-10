-- Deploy instrument_data_attribute_permission
-- requires: instrument_data_attribute

BEGIN;

REVOKE ALL ON TABLE instrument.data_attribute FROM PUBLIC;
REVOKE ALL ON TABLE instrument.data_attribute FROM genome;
GRANT ALL ON TABLE instrument.data_attribute TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE instrument.data_attribute TO "gms-user";

COMMIT;
