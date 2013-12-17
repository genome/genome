-- Revert web_talk

BEGIN;

DROP TABLE IF EXISTS web.task;

COMMIT;
