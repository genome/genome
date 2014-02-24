-- Verify config.analysismenu_item.description

BEGIN;

ALTER TABLE config.analysismenu_item DROP COLUMN description;

ROLLBACK;
