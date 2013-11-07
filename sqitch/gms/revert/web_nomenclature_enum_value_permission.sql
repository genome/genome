-- Revert web_nomenclature_enum_value_permission

BEGIN;

REVOKE ALL ON TABLE web.nomenclature_enum_value FROM "gms-user";

COMMIT;
