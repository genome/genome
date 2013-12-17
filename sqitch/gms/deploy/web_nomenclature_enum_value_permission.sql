-- Deploy web_nomenclature_enum_value_permission
-- requires: web_nomenclature_enum_value

BEGIN;

REVOKE ALL ON TABLE web.nomenclature_enum_value FROM PUBLIC;
REVOKE ALL ON TABLE web.nomenclature_enum_value FROM genome;
GRANT ALL ON TABLE web.nomenclature_enum_value TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE web.nomenclature_enum_value TO "gms-user";

COMMIT;
