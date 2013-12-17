-- Revert web_nomenclature_field_permission

BEGIN;

REVOKE ALL ON TABLE web.nomenclature_field FROM "gms-user";

COMMIT;
