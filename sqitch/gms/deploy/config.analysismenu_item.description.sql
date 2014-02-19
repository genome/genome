-- Deploy config.analysismenu_item.description
-- requires: config_analysismenu_item

BEGIN;

ALTER TABLE config.analysismenu_item ADD COLUMN description TEXT;

COMMIT;
